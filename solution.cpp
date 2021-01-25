#include <iostream> 
#include "functions.h"
#include <math.h>
#include <sys/time.h>
using namespace std;

int solution( int n, int m, double *a, double *b, double *x, double *v_1, double *v_2)
{
	int k , l, last, min, s = 0, j = 0, i = 0, a_v, a_h, t = 0, s1, q, last1, d, s_n;
	double c00 = 0.0, c01 = 0.0, c02 = 0.0, c10 = 0.0, c11 = 0.0, c12 = 0.0, c20 = 0.0, c21 = 0.0, c22 = 0.0;
	double *p_a;
	double sum = 0.0;
	double eps = 1.0;
	double divisor = 0.0;
	double norm = 0.0;
	int r = 0, z = 0, y = 0, r_m = 0, z_m = 0, m_new = 0;

	for(i = 0; i < m * m; i ++)
		v_2[i] = 0.0;

	eps = 1.0;
	while( 1.0 + eps > 1.0 )
	{
		eps = eps / 2.0;
	}

	// column norm 

	for(i = 0; i < n; i ++) norm += fabs(a[i*n]);

	for(j = 1; j < n; j ++)
	{
		for(i = 0; i < n; i ++) sum += fabs(a[i*n + j]);

		if( sum > norm ) norm = sum;

		sum = 0.0;
	}

	sum = 0.0;

	//struct timeval tv1;
	//struct timeval tv2;
 	//struct timezone tz;
  	//gettimeofday(&tv1, &tz);

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

	if(m % 3 == 0)
	{
		m_new = m;

		if(l == 0)
		{
			a_v = m;
			a_h = m;
	
			for(s = 0; s < last; s ++)
			{
				// BEGIN making A_ss E
	
				s_n = s*n;
		
				p_a = a + (s_n+s)*m;
		
				last1 = m;
		
				
					//return_i_j_matrix_block(s, s, n, m, k, l, a, v_1);// put A_ss into V_1
		
					//inverse_matrix(m, v_1, v_2); // V_2 = ( A_ss )^(-1)
		
				divisor = 0.0;
		
				for(q = 0; q < m; q ++ )
				{
					for( d = 0; d < m; d ++ ) 
					{
						if( d == q) v_1[q*m + d] = 1.0;
						else v_1[q*m + d] = 0.0;
					}
				}
		
				for(d = 0; d < m; d++ )  // step number s
				{	
					//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
					if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
					{
						printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
						printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
						return 1;
					}
		
					// futher a_ss != 0
					divisor = p_a[d*n + d];
		
					for(t = 0; t < m; t++) // made a_ss equal to 1
					{
						p_a[d*n + t] = p_a[d*n + t] / divisor; 
		
						v_1[d*m + t] = v_1[d*m + t] / divisor; 
					}
		
					for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
					for(s1 = d + 1; s1 < m; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
				}
		
				//t = 0;
		
				//matrix_multiplication(m, m, 1, v_2, v_1, v_3);// ( A_ss )^(-1) * B_s
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);;
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
					//put_back_i_vector_part(s, m, k, l, b, v_3); // V_2, V_3 and V_1 are free
		
					/*if( l != 0)
					{
						if( s == k) last1 = l;
						else last1 = m;
					}
					else last1 = m;*/
		
				for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
		
				// END making A_ss E
				// BEGIN ( A_ss )^(-1) * A_sj
		
				for( j = s + 1; j < last; j ++ )
				{
					p_a = a + (s_n+j)*m;
					
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_1);// put A_ij into V_1
		
					//matrix_multiplication(m, m, m, v_2, v_1, v_3);// V_3 = ( A_ss )^(-1) * A_sj
		
					/*for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < m; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]);
		
							v_2[s1*m + t] = sum;
						}
					}*/


					//print_matrix(9, 9, 9, v_2);


					for(r = 0; r < m; r += 3)
					{
						//printf("HEREEEEEEEE\n");
						r_m = r * m;
				
						for(y = 0; y < m; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
				
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
				
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
				
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
								c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
							}
				
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
							
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
				
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
							v_2[r_m + 2*m + y + 2] = c22;
				
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
				
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
				
							c20 = 0.0;
							c21 = 0.0;
							c22 = 0.0;
						}
					}

					//print_matrix(9, 9, 9, v_2);










					//t = 0;
		
					//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3);// V_3 -> A_sj
		
					/*if( s == k) a_v = l;
					else a_v = m;
					if( j == k) a_h = l;
					else a_h = m;*/
		
					p_a = a + (s_n+j)*m;
	
					t = 0;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							p_a[ s1*n + q ] = v_2[t];
							t++;
						}
					}
					//t = 0;
					
				}
		
				// END ( A_ss )^(-1) * A_sj
				// BEGIN do row difference 
		
				for( i = s + 1; i < last; i ++ )
				{
		
					//return_i_j_matrix_block(i, s, n, m, k, l, a, v_1);// put A_is into V_1
		
					t = 0;
		
					p_a = a + (i*n+s)*m;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							v_1[t] = p_a[ s1*n + q ];
							p_a[ s1*n + q ] = 0.0;
							t++;
						}
					}
					//t = 0;
		
					//matrix_multiplication(m, m, 1, v_1, v_2, v_2); 
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
		
						//matrix_difference(m, 1, v_2, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
		
					for( j = s + 1; j < last; j ++ )
					{
						//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (s_n+j)*m;
		
						//matrix_multiplication(m, m, m, v_1, v_2, v_3);
		
						/*for(s1 = 0; s1 < m; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
		
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
		
								v_2[s1*m + t] = sum;
							}
						}*/
						//t = 0;

						for(r = 0; r < m; r += 3)
						{
							//printf("HEREEEEEEEE\n");
							r_m = r * m;
					
							for(y = 0; y < m; y += 3)
							{
								for(z = 0; z < m; z ++)
								{
									z_m = z * n;
					
									c00 += v_1[r_m + z] * p_a[z_m + y];
									c01 += v_1[r_m + z] * p_a[z_m + y + 1];
									c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
									c10 += v_1[r_m + m + z] * p_a[z_m + y];
									c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
									c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
									c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
									c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
									c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
								}
					
								v_2[r_m + y]     = c00;
								v_2[r_m + y + 1] = c01;
								v_2[r_m + y + 2] = c02;
								
								v_2[r_m + m + y]     = c10;
								v_2[r_m + m + y + 1] = c11;
								v_2[r_m + m + y + 2] = c12;
					
								v_2[r_m + 2*m + y]     = c20;
								v_2[r_m + 2*m + y + 1] = c21;
								v_2[r_m + 2*m + y + 2] = c22;
					
								c00 = 0.0;
								c01 = 0.0;
								c02 = 0.0;
					
								c10 = 0.0;
								c11 = 0.0;
								c12 = 0.0;
					
								c20 = 0.0;
								c21 = 0.0;
								c22 = 0.0;
							}
						}





















		
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (i*n+j)*m;
		
						//matrix_difference(m, m, v_2, v_3);
		
						for(s1 = 0; s1 < m; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					}
				}
				// END do row difference
			}
			
			// BEGIN of X_i
		
			for( i = 0; i < min; i ++)	x[(last - 1)*m + i] = b[(last - 1)*m + i];
		
			for( i = last - 2; i >= 0; i --)
			{
		
				//return_i_vector_part(last - 1, m, k, l, b, v_1);
		
				p_a = a + (i*n+(last - 1))*m;
		
				// matrix_multiplication(m, m, 1, v_2, v_1, v_3);
				
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[(last - 1)*m + q*1 + t]);//(v_2[s1*m + q]) * (v_1[q*1 + t]);
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
				//matrix_difference(m, 1, v_1, v_3);
		
				for(s1 = 0; s1 < m; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
				//t = 0;
		
				for( j = last - 2; j > i ; j --)
				{ 
		
					//return_i_vector_part(j, m, k, l, x, v_1);// put X_i into V_2
		
					p_a = a + (i*n+j)*m;
		
					//matrix_multiplication(m, m, 1, v_2, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[j*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
					//t = 0;
		
					//matrix_difference(m, 1, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
				}
		
				for( j = 0; j < m; j ++)	x[i*m + j] = b[i*m + j];
		
			}
		}
	
		else // l != 0
		{
			a_v = m;
			a_h = m;
	
			for(s = 0; s < last - 1; s ++)
			{
				// BEGIN making A_ss E
	
				s_n = s * n;
		
				p_a = a + (s_n+s)*m;
		
				last1 = m;
		
				
					//return_i_j_matrix_block(s, s, n, m, k, l, a, v_1);// put A_ss into V_1
		
					//inverse_matrix(m, v_1, v_2); // V_2 = ( A_ss )^(-1)
		
				for(q = 0; q < m; q ++ )
				{
					for( d = 0; d < m; d ++ ) 
					{
						if( d == q) v_1[q*m + d] = 1.0;
						else v_1[q*m + d] = 0.0;
					}
				}
		
				for(d = 0; d < m; d++ )  // step number s
				{	
					//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
					if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
					{
						printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
						printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
						return 1;
					}
		
					// futher a_ss != 0
					divisor = p_a[d*n + d];
		
					for(t = 0; t < m; t++) // made a_ss equal to 1
					{
						p_a[d*n + t] = p_a[d*n + t] / divisor; 
		
						v_1[d*m + t] = v_1[d*m + t] / divisor; 
					}
		
					for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
					for(s1 = d + 1; s1 < m; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
				}
		
				//t = 0;
		
				//matrix_multiplication(m, m, 1, v_2, v_1, v_3);// ( A_ss )^(-1) * B_s
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);;
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
					//put_back_i_vector_part(s, m, k, l, b, v_3); // V_2, V_3 and V_1 are free
		
					/*if( l != 0)
					{
						if( s == k) last1 = l;
						else last1 = m;
					}
					else last1 = m;*/
		
				for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
		
				// END making A_ss E
				// BEGIN ( A_ss )^(-1) * A_sj
		
				for( j = s + 1; j < last - 1; j ++ )
				{
					p_a = a + (s_n+j)*m;
					
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_1);// put A_ij into V_1
		
					//matrix_multiplication(m, m, m, v_2, v_1, v_3);// V_3 = ( A_ss )^(-1) * A_sj
		
					/*for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < m; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]);
		
							v_2[s1*m + t] = sum;
						}
					}*/






					for(r = 0; r < m; r += 3)
					{
						//printf("HEREEEEEEEE\n");
						r_m = r * m;
				
						for(y = 0; y < m; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
				
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
				
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
				
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
								c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
							}
				
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
							
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
				
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
							v_2[r_m + 2*m + y + 2] = c22;
				
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
				
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
				
							c20 = 0.0;
							c21 = 0.0;
							c22 = 0.0;
						}
					}
					//t = 0;
		
					//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3);// V_3 -> A_sj
		
					/*if( s == k) a_v = l;
					else a_v = m;
					if( j == k) a_h = l;
					else a_h = m;*/
		
					//p_a = a + (s*n+j)*m;
					t = 0;
		
					for( s1 = 0; s1 < m; s1++)
					{
						for(q = 0; q < m; q++) 
						{
							p_a[ s1*n + q ] = v_2[t];
							t++;
						}
					}
					//t = 0;
					
				}
	
				// j == last - 1
				j = last - 1;
	
				p_a = a + (s_n + j)*m;
	
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < min; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]); //(v_1[q*min + t]);
	
						v_2[s1*min + t] = sum;
					}
				}
				//t = 0;
	
				//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3); // V_3 -> A_sj
	
				/*if( s == k) a_v = l;
				else a_v = m;
				if( j == k) a_h = l;
				else a_h = m;*/
	
				p_a = a + (s_n + j)*m;
	
				t = 0;
	
				for( s1 = 0; s1 < m; s1++)
				{
					for(q = 0; q < l; q++) 
					{
						p_a[ s1*n + q ] = v_2[t];
						t++;
					}
				}
				//t = 0;
	
	
	
	
	
		
				// END ( A_ss )^(-1) * A_sj
				// BEGIN do row difference 
	
				a_h = m;
		
				for( i = s + 1; i < last - 1; i ++ )
				{
	
					//print_matrix(n, n, 7, a);
		
					//return_i_j_matrix_block(i, s, n, m, k, l, a, v_1);// put A_is into V_1
		
					t = 0;
		
					p_a = a + (i*n+s)*m;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							v_1[t] = p_a[ s1*n + q ];
							p_a[ s1*n + q ] = 0.0;
							t++;
						}
					}
					//t = 0;
		
					//matrix_multiplication(m, m, 1, v_1, v_2, v_2); 
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
		
						//matrix_difference(m, 1, v_2, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
		
					for( j = s + 1; j < last - 1; j ++ )
					{
						//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (s_n + j)*m;
		
						//matrix_multiplication(m, m, m, v_1, v_2, v_3);
		
						/*for(s1 = 0; s1 < m; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
		
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
		
								v_2[s1*m + t] = sum;
							}
						}*/


						for(r = 0; r < m; r += 3)
						{
							//printf("HEREEEEEEEE\n");
							r_m = r * m;
					
							for(y = 0; y < m; y += 3)
							{
								for(z = 0; z < m; z ++)
								{
									z_m = z * n;
					
									c00 += v_1[r_m + z] * p_a[z_m + y];
									c01 += v_1[r_m + z] * p_a[z_m + y + 1];
									c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
									c10 += v_1[r_m + m + z] * p_a[z_m + y];
									c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
									c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
									c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
									c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
									c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
								}
					
								v_2[r_m + y]     = c00;
								v_2[r_m + y + 1] = c01;
								v_2[r_m + y + 2] = c02;
								
								v_2[r_m + m + y]     = c10;
								v_2[r_m + m + y + 1] = c11;
								v_2[r_m + m + y + 2] = c12;
					
								v_2[r_m + 2*m + y]     = c20;
								v_2[r_m + 2*m + y + 1] = c21;
								v_2[r_m + 2*m + y + 2] = c22;
					
								c00 = 0.0;
								c01 = 0.0;
								c02 = 0.0;
					
								c10 = 0.0;
								c11 = 0.0;
								c12 = 0.0;
					
								c20 = 0.0;
								c21 = 0.0;
								c22 = 0.0;
							}
						}








						//t = 0;
		
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (i*n+j)*m;
		
						//matrix_difference(m, m, v_2, v_3);
		
						for(s1 = 0; s1 < m; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					}
	
					//last column for not last row
	
					j = last - 1;
	
					p_a = a + (s_n + j)*m;
	
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < min; t ++ )
						{
							sum = 0.0;
	
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
							v_2[s1*min + t] = sum;
						}
					}
					//t = 0;
	
					//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
					p_a = a + (i*n+j)*m;
	
					//matrix_difference(m, min, v_2, v_3);
	
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*min + t];
	
					//t = 0;
				}
	
	
				//row difference for last row
				i = last - 1;
	
				p_a = a + (i*n+s)*m;
	
				t = 0;
	
				for( s1 = 0; s1 < l; s1++)
				{
					for(q = 0; q < m; q++) 
					{
						v_1[t] = p_a[ s1*n + q ];
						p_a[ s1*n + q ] = 0.0;
						t++;
					}
				}
				//t = 0;
	
				for(s1 = 0; s1 < min; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
	
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
	
				//matrix_difference(min, 1, v_1, v_3);// V_2 = V_2 - V_3
	
				for(s1 = 0; s1 < min; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
	
				//t = 0;
	
				for( j = s + 1; j < last - 1; j ++ )
				{
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
					//p_a = a + (s*n+j)*m;
	
						//matrix_multiplication(min, m, m, v_1, v_2, v_3);
	
						p_a = a + (s_n + j)*m;
	
						for(s1 = 0; s1 < min; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
	
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
								v_2[s1*m + t] = sum;
							}
						}
						//t = 0;
	
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
						p_a = a + (i*n+j)*m;
	
						//matrix_difference(min, m, v_2, v_3);
	
						for(s1 = 0; s1 < min; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					
				}
	
				j = last - 1;
					
						//matrix_multiplication(min, m, min, v_1, v_2, v_3);
				p_a = a + (s_n + j)*m;
	
				for(s1 = 0; s1 < min; s1 ++)
				{
					for(t = 0; t < min; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
						v_2[s1*min + t] = sum;
					}
				}
	
				p_a = a + (i*n+j)*m;
	
				//matrix_difference(min, min, v_2, v_3);
	
				for(s1 = 0; s1 < min; s1 ++)
					for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*min + t];
	
				//t = 0;
				// END do row difference
			}
	
			// inverse for last row
			s = last - 1;
	
			s_n = s * n;
	
			p_a = a + (s_n +s )*m;
	
			last1 = l;
	
			for(q = 0; q < min; q ++ )
			{
				for( d = 0; d < min; d ++ ) 
				{
					if( d == q) v_1[q*min + d] = 1.0;
					else v_1[q*min + d] = 0.0;
				}
			}
	
			for(d = 0; d < min; d++ )  // step number s
			{
				//print_matrix(n, n, r, a);
				//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
				if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
				{
					printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
					printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
					return 1;
				}
	
				// futher a_ss != 0
				divisor = p_a[d*n + d];
				//printf("\ndivisor is %10.3e\n", divisor);
	
				for(t = 0; t < min; t++) // made a_ss equal to 1
				{
					p_a[d*n + t] = p_a[d*n + t] / divisor; 
	
					v_1[d*min + t] = v_1[d*min + t] / divisor; 
				}
	
				//print_matrix(n, n, r, v_1);
	
				for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
				{
					//if( s1 != d) 
					//{
						if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
						{
							divisor = p_a[s1*n + d];
	
							//printf("\ncycle divisor is %10.3e\n", divisor);
	
							for( t = 0; t < min; t++ )	
							{
								//printf("\n%10.10Le = %10.20Le - %10.20Le * %10.20Le\n", v_1[s1*min + t] - divisor * v_1[d*min + t], v_1[s1*min + t], divisor, v_1[d*min + t]);
	
								p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
	
								v_1[s1*min + t] = v_1[s1*min + t] - divisor * v_1[d*min + t];
							}
	
						}
					//}
				}
				for(s1 = d + 1; s1 < min; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
				{
					//if( s1 != d) 
					//{
						if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
						{
							divisor = p_a[s1*n + d];
	
							//printf("\ncycle divisor is %10.3e\n", divisor);
	
							for( t = 0; t < min; t++ )	
							{
								//printf("\n%10.10Le = %10.20Le - %10.20Le * %10.20Le\n", v_1[s1*min + t] - divisor * v_1[d*min + t], v_1[s1*min + t], divisor, v_1[d*min + t]);
	
								p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
	
								v_1[s1*min + t] = v_1[s1*min + t] - divisor * v_1[d*min + t];
							}
	
						}
					//}
				}
			}
	
			//t = 0;
	
			for(s1 = 0; s1 < min; s1 ++)
			{
				for(t = 0; t < 1; t ++ )
				{
					sum = 0.0;
	
					for(q = 0; q < min; q ++ ) sum += (v_1[s1*min + q]) * (b[s*m + q*1 + t]);
	
					v_2[s1*1 + t] = sum;
				}
		    }
			//t = 0;
	
			//put_back_i_vector_part(s, m, k, l, b, v_3); // V_3 and V_1 are free
	
			/*if( l != 0)
			{
				if( s == k) last1 = l;
				else last1 = m;
			}
			else last1 = m;*/
	
			for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
	
			//print_matrix(n, n, 7, a);
	
			//for(i = 0; i < 7; i++) printf(" %10.3e\n", b[i]);
	
			// BEGIN of X_i
		
			for( i = 0; i < min; i ++)	x[(last - 1)*m + i] = b[(last - 1)*m + i];
		
			for( i = last - 2; i >= 0; i --)
			{
		
				//return_i_vector_part(last - 1, m, k, l, b, v_1);
		
				p_a = a + (i*n+(last - 1))*m;
		
				// matrix_multiplication(m, m, 1, v_2, v_1, v_3);
				
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < min; q ++ ) sum += (p_a[s1*n + q]) * (b[(last - 1)*m + q*1 + t]);//(v_2[s1*m + q]) * (v_1[q*1 + t]);
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
				//matrix_difference(m, 1, v_1, v_3);
		
				for(s1 = 0; s1 < m; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
				//t = 0;
		
				for( j = last - 2; j > i ; j --)
				{ 
		
					//return_i_vector_part(j, m, k, l, x, v_1);// put X_i into V_2
		
					p_a = a + (i*n+j)*m;
		
					//matrix_multiplication(m, m, 1, v_2, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[j*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
					//t = 0;
		
					//matrix_difference(m, 1, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
				}
		
				for( j = 0; j < m; j ++)	x[i*m + j] = b[i*m + j];
		
			}
		}	
	}








/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







	if(m % 3 == 1)
	{
		//printf("HEREEEEEEEE\n");
		m_new = m - 1;

		if(l == 0)
		{
			a_v = m;
			a_h = m;
	
			for(s = 0; s < last; s ++)
			{
				// BEGIN making A_ss E
	
				s_n = s*n;
		
				p_a = a + (s_n+s)*m;
		
				last1 = m;
		
				
					//return_i_j_matrix_block(s, s, n, m, k, l, a, v_1);// put A_ss into V_1
		
					//inverse_matrix(m, v_1, v_2); // V_2 = ( A_ss )^(-1)
		
				divisor = 0.0;
		
				for(q = 0; q < m; q ++ )
				{
					for( d = 0; d < m; d ++ ) 
					{
						if( d == q) v_1[q*m + d] = 1.0;
						else v_1[q*m + d] = 0.0;
					}
				}
		
				for(d = 0; d < m; d++ )  // step number s
				{	
					//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
					if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
					{
						printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
						printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
						return 1;
					}
		
					// futher a_ss != 0
					divisor = p_a[d*n + d];
		
					for(t = 0; t < m; t++) // made a_ss equal to 1
					{
						p_a[d*n + t] = p_a[d*n + t] / divisor; 
		
						v_1[d*m + t] = v_1[d*m + t] / divisor; 
					}
		
					for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
					for(s1 = d + 1; s1 < m; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
				}
		
				//t = 0;
		
				//matrix_multiplication(m, m, 1, v_2, v_1, v_3);// ( A_ss )^(-1) * B_s
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);;
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
					//put_back_i_vector_part(s, m, k, l, b, v_3); // V_2, V_3 and V_1 are free
		
					/*if( l != 0)
					{
						if( s == k) last1 = l;
						else last1 = m;
					}
					else last1 = m;*/
		
				for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
		
				// END making A_ss E
				// BEGIN ( A_ss )^(-1) * A_sj
		
				for( j = s + 1; j < last; j ++ )
				{
					p_a = a + (s_n+j)*m;
					
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_1);// put A_ij into V_1
		
					//matrix_multiplication(m, m, m, v_2, v_1, v_3);// V_3 = ( A_ss )^(-1) * A_sj
		
					/*for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < m; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]);
		
							v_2[s1*m + t] = sum;
						}
					}*/


					//print_matrix(9, 9, 9, v_2);

					for(r = 0; r < m_new; r += 3)
					{
						r_m = r * m;
				
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
				
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
				
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
				
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
								c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
							}
				
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
							
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
				
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
							v_2[r_m + 2*m + y + 2] = c22;
				
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
				
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
				
							c20 = 0.0;
							c21 = 0.0;
							c22 = 0.0;
						}
				
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
				
							c00 += v_1[r_m + z] * p_a[z_m + y];
				
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
				
							c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
						}
				
						v_2[r_m + y]     = c00;
						
						v_2[r_m + m + y]     = c10;
						
						v_2[r_m + 2*m + y]     = c20;
						
						c00 = 0.0;
						
						c10 = 0.0;
						
						c20 = 0.0;
					}
				
					//for(r = m_new; r < m; r ++)
					//{
					r_m = m_new * m;
					
					for(y = 0; y < m_new; y += 3)
					{
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
					
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
							c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
							//c10 += v_1[r_m + m + z] * p_a[z_m + y];
							//c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
							//c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
						}
					
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
						v_2[r_m + y + 2] = c02;
								
						//v_2[r_m + m + y]     = c10;
						//v_2[r_m + m + y + 1] = c11;
						//v_2[r_m + m + y + 2] = c12;
					
						c00 = 0.0;
						c01 = 0.0;
						c02 = 0.0;
					}
						//}
					for(z = 0; z < m; z ++)
					{
						z_m = z * n;
					
						c00 += v_1[r_m + z] * p_a[z_m + y];
					
						c10 += v_1[r_m + m + z] * p_a[z_m + y];
					}
					
					v_2[r_m + y] = c00;
					
					//v_2[r_m + m + y] = c10;
					
					c00 = 0.0;
					
					c10 = 0.0;




					//t = 0;
		
					//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3);// V_3 -> A_sj
		
					/*if( s == k) a_v = l;
					else a_v = m;
					if( j == k) a_h = l;
					else a_h = m;*/
		
					p_a = a + (s_n+j)*m;
	
					t = 0;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							p_a[ s1*n + q ] = v_2[t];
							t++;
						}
					}
					//t = 0;
					
				}
		
				// END ( A_ss )^(-1) * A_sj
				// BEGIN do row difference 
		
				for( i = s + 1; i < last; i ++ )
				{
		
					//return_i_j_matrix_block(i, s, n, m, k, l, a, v_1);// put A_is into V_1
		
					t = 0;
		
					p_a = a + (i*n+s)*m;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							v_1[t] = p_a[ s1*n + q ];
							p_a[ s1*n + q ] = 0.0;
							t++;
						}
					}
					//t = 0;
		
					//matrix_multiplication(m, m, 1, v_1, v_2, v_2); 
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
		
						//matrix_difference(m, 1, v_2, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
		
					for( j = s + 1; j < last; j ++ )
					{
						//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (s_n+j)*m;
		
						//matrix_multiplication(m, m, m, v_1, v_2, v_3);
		
						/*for(s1 = 0; s1 < m; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
		
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
		
								v_2[s1*m + t] = sum;
							}
						}*/
						//t = 0;

						for(r = 0; r < m_new; r += 3)
						{
							r_m = r * m;
					
							for(y = 0; y < m_new; y += 3)
							{
								for(z = 0; z < m; z ++)
								{
									z_m = z * n;
					
									c00 += v_1[r_m + z] * p_a[z_m + y];
									c01 += v_1[r_m + z] * p_a[z_m + y + 1];
									c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
									c10 += v_1[r_m + m + z] * p_a[z_m + y];
									c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
									c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
									c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
									c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
									c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
								}
					
								v_2[r_m + y]     = c00;
								v_2[r_m + y + 1] = c01;
								v_2[r_m + y + 2] = c02;
								
								v_2[r_m + m + y]     = c10;
								v_2[r_m + m + y + 1] = c11;
								v_2[r_m + m + y + 2] = c12;
					
								v_2[r_m + 2*m + y]     = c20;
								v_2[r_m + 2*m + y + 1] = c21;
								v_2[r_m + 2*m + y + 2] = c22;
					
								c00 = 0.0;
								c01 = 0.0;
								c02 = 0.0;
					
								c10 = 0.0;
								c11 = 0.0;
								c12 = 0.0;
					
								c20 = 0.0;
								c21 = 0.0;
								c22 = 0.0;
							}
					
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
					
								c00 += v_1[r_m + z] * p_a[z_m + y];
					
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
					
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
							}
					
							v_2[r_m + y]     = c00;
							
							v_2[r_m + m + y]     = c10;
							
							v_2[r_m + 2*m + y]     = c20;
							
							c00 = 0.0;
							
							c10 = 0.0;
							
							c20 = 0.0;
						}
					
						//for(r = m_new; r < m; r ++)
						//{
						r_m = m_new * m;
						
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
						
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
						
								//c10 += v_1[r_m + m + z] * p_a[z_m + y];
								//c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								//c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
						
							}
						
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
									
							//v_2[r_m + m + y]     = c10;
							//v_2[r_m + m + y + 1] = c11;
							//v_2[r_m + m + y + 2] = c12;
						
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
						}
							//}
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
						
							c00 += v_1[r_m + z] * p_a[z_m + y];
						
							//c10 += v_1[r_m + m + z] * p_a[z_m + y];
						}
						
						v_2[r_m + y]     = c00;
						
						//v_2[r_m + m + y]     = c10;
						
						c00 = 0.0;
						
						//c10 = 0.0;





















		
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (i*n+j)*m;
		
						//matrix_difference(m, m, v_2, v_3);
		
						for(s1 = 0; s1 < m; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					}
				}
				// END do row difference
			}
			
			// BEGIN of X_i
		
			for( i = 0; i < min; i ++)	x[(last - 1)*m + i] = b[(last - 1)*m + i];
		
			for( i = last - 2; i >= 0; i --)
			{
		
				//return_i_vector_part(last - 1, m, k, l, b, v_1);
		
				p_a = a + (i*n+(last - 1))*m;
		
				// matrix_multiplication(m, m, 1, v_2, v_1, v_3);
				
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[(last - 1)*m + q*1 + t]);//(v_2[s1*m + q]) * (v_1[q*1 + t]);
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
				//matrix_difference(m, 1, v_1, v_3);
		
				for(s1 = 0; s1 < m; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
				//t = 0;
		
				for( j = last - 2; j > i ; j --)
				{ 
		
					//return_i_vector_part(j, m, k, l, x, v_1);// put X_i into V_2
		
					p_a = a + (i*n+j)*m;
		
					//matrix_multiplication(m, m, 1, v_2, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[j*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
					//t = 0;
		
					//matrix_difference(m, 1, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
				}
		
				for( j = 0; j < m; j ++)	x[i*m + j] = b[i*m + j];
		
			}
		}
	
		else // l != 0
		{
			a_v = m;
			a_h = m;
	
			for(s = 0; s < last - 1; s ++)
			{
				// BEGIN making A_ss E
	
				s_n = s * n;
		
				p_a = a + (s_n+s)*m;
		
				last1 = m;
		
				
					//return_i_j_matrix_block(s, s, n, m, k, l, a, v_1);// put A_ss into V_1
		
					//inverse_matrix(m, v_1, v_2); // V_2 = ( A_ss )^(-1)
		
				for(q = 0; q < m; q ++ )
				{
					for( d = 0; d < m; d ++ ) 
					{
						if( d == q) v_1[q*m + d] = 1.0;
						else v_1[q*m + d] = 0.0;
					}
				}
		
				for(d = 0; d < m; d++ )  // step number s
				{	
					//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
					if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
					{
						printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
						printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
						return 1;
					}
		
					// futher a_ss != 0
					divisor = p_a[d*n + d];
		
					for(t = 0; t < m; t++) // made a_ss equal to 1
					{
						p_a[d*n + t] = p_a[d*n + t] / divisor; 
		
						v_1[d*m + t] = v_1[d*m + t] / divisor; 
					}
		
					for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
					for(s1 = d + 1; s1 < m; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
				}
		
				//t = 0;
		
				//matrix_multiplication(m, m, 1, v_2, v_1, v_3);// ( A_ss )^(-1) * B_s
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);;
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
					//put_back_i_vector_part(s, m, k, l, b, v_3); // V_2, V_3 and V_1 are free
		
					/*if( l != 0)
					{
						if( s == k) last1 = l;
						else last1 = m;
					}
					else last1 = m;*/
		
				for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
		
				// END making A_ss E
				// BEGIN ( A_ss )^(-1) * A_sj
		
				for( j = s + 1; j < last - 1; j ++ )
				{
					p_a = a + (s_n+j)*m;
					
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_1);// put A_ij into V_1
		
					//matrix_multiplication(m, m, m, v_2, v_1, v_3);// V_3 = ( A_ss )^(-1) * A_sj
		
					/*for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < m; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]);
		
							v_2[s1*m + t] = sum;
						}
					}*/

					for(r = 0; r < m_new; r += 3)
					{
						r_m = r * m;
				
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
				
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
				
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
				
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
								c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
							}
				
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
							
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
				
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
							v_2[r_m + 2*m + y + 2] = c22;
				
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
				
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
				
							c20 = 0.0;
							c21 = 0.0;
							c22 = 0.0;
						}
				
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
				
							c00 += v_1[r_m + z] * p_a[z_m + y];
				
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
				
							c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
						}
				
						v_2[r_m + y]     = c00;
						
						v_2[r_m + m + y]     = c10;
						
						v_2[r_m + 2*m + y]     = c20;
						
						c00 = 0.0;
						
						c10 = 0.0;
						
						c20 = 0.0;
					}
				
					//for(r = m_new; r < m; r ++)
					//{
					r_m = m_new * m;
					
					for(y = 0; y < m_new; y += 3)
					{
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
					
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
							c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
							//c10 += v_1[r_m + m + z] * p_a[z_m + y];
							//c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
							//c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
						}
					
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
						v_2[r_m + y + 2] = c02;
								
						//v_2[r_m + m + y]     = c10;
						//v_2[r_m + m + y + 1] = c11;
						//v_2[r_m + m + y + 2] = c12;
					
						c00 = 0.0;
						c01 = 0.0;
						c02 = 0.0;
					}
						//}
					for(z = 0; z < m; z ++)
					{
						z_m = z * n;
					
						c00 += v_1[r_m + z] * p_a[z_m + y];
					
						//c10 += v_1[r_m + m + z] * p_a[z_m + y];
					}
					
					v_2[r_m + y]     = c00;
					
					//v_2[r_m + m + y]     = c10;
					
					c00 = 0.0;
					
					//c10 = 0.0;












					//t = 0;
		
					//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3);// V_3 -> A_sj
		
					/*if( s == k) a_v = l;
					else a_v = m;
					if( j == k) a_h = l;
					else a_h = m;*/
		
					//p_a = a + (s*n+j)*m;
					t = 0;
		
					for( s1 = 0; s1 < m; s1++)
					{
						for(q = 0; q < m; q++) 
						{
							p_a[ s1*n + q ] = v_2[t];
							t++;
						}
					}
					//t = 0;
					
				}
	
				// j == last - 1
				j = last - 1;
	
				p_a = a + (s_n + j)*m;
	
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < min; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]); //(v_1[q*min + t]);
	
						v_2[s1*min + t] = sum;
					}
				}
				//t = 0;
	
				//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3); // V_3 -> A_sj
	
				/*if( s == k) a_v = l;
				else a_v = m;
				if( j == k) a_h = l;
				else a_h = m;*/
	
				p_a = a + (s_n + j)*m;
	
				t = 0;
	
				for( s1 = 0; s1 < m; s1++)
				{
					for(q = 0; q < l; q++) 
					{
						p_a[ s1*n + q ] = v_2[t];
						t++;
					}
				}
				//t = 0;
	
	
	
	
	
		
				// END ( A_ss )^(-1) * A_sj
				// BEGIN do row difference 
	
				a_h = m;
		
				for( i = s + 1; i < last - 1; i ++ )
				{
	
					//print_matrix(n, n, 7, a);
		
					//return_i_j_matrix_block(i, s, n, m, k, l, a, v_1);// put A_is into V_1
		
					t = 0;
		
					p_a = a + (i*n+s)*m;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							v_1[t] = p_a[ s1*n + q ];
							p_a[ s1*n + q ] = 0.0;
							t++;
						}
					}
					//t = 0;
		
					//matrix_multiplication(m, m, 1, v_1, v_2, v_2); 
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
		
						//matrix_difference(m, 1, v_2, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
		
					for( j = s + 1; j < last - 1; j ++ )
					{
						//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (s_n + j)*m;
		
						//matrix_multiplication(m, m, m, v_1, v_2, v_3);
		
						/*for(s1 = 0; s1 < m; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
		
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
		
								v_2[s1*m + t] = sum;
							}
						}*/


						for(r = 0; r < m_new; r += 3)
						{
							r_m = r * m;
					
							for(y = 0; y < m_new; y += 3)
							{
								for(z = 0; z < m; z ++)
								{
									z_m = z * n;
					
									c00 += v_1[r_m + z] * p_a[z_m + y];
									c01 += v_1[r_m + z] * p_a[z_m + y + 1];
									c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
									c10 += v_1[r_m + m + z] * p_a[z_m + y];
									c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
									c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
									c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
									c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
									c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
								}
					
								v_2[r_m + y]     = c00;
								v_2[r_m + y + 1] = c01;
								v_2[r_m + y + 2] = c02;
								
								v_2[r_m + m + y]     = c10;
								v_2[r_m + m + y + 1] = c11;
								v_2[r_m + m + y + 2] = c12;
					
								v_2[r_m + 2*m + y]     = c20;
								v_2[r_m + 2*m + y + 1] = c21;
								v_2[r_m + 2*m + y + 2] = c22;
					
								c00 = 0.0;
								c01 = 0.0;
								c02 = 0.0;
					
								c10 = 0.0;
								c11 = 0.0;
								c12 = 0.0;
					
								c20 = 0.0;
								c21 = 0.0;
								c22 = 0.0;
							}
					
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
					
								c00 += v_1[r_m + z] * p_a[z_m + y];
					
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
					
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
							}
					
							v_2[r_m + y]     = c00;
							
							v_2[r_m + m + y]     = c10;
							
							v_2[r_m + 2*m + y]     = c20;
							
							c00 = 0.0;
							
							c10 = 0.0;
							
							c20 = 0.0;
						}
					
						//for(r = m_new; r < m; r ++)
						//{
						r_m = m_new * m;
						
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
						
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
						
								//c10 += v_1[r_m + m + z] * p_a[z_m + y];
								//c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								//c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
						
							}
						
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
									
							//v_2[r_m + m + y]     = c10;
							//v_2[r_m + m + y + 1] = c11;
							//v_2[r_m + m + y + 2] = c12;
						
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
						}
							//}
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
						
							c00 += v_1[r_m + z] * p_a[z_m + y];
						
							//c10 += v_1[r_m + m + z] * p_a[z_m + y];
						}
						
						v_2[r_m + y]     = c00;
						
						//v_2[r_m + m + y]     = c10;
						
						c00 = 0.0;
						
						//c10 = 0.0;








						//t = 0;
		
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (i*n+j)*m;
		
						//matrix_difference(m, m, v_2, v_3);
		
						for(s1 = 0; s1 < m; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					}
	
					//last column for not last row
	
					j = last - 1;
	
					p_a = a + (s_n + j)*m;
	
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < min; t ++ )
						{
							sum = 0.0;
	
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
							v_2[s1*min + t] = sum;
						}
					}
					//t = 0;
	
					//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
					p_a = a + (i*n+j)*m;
	
					//matrix_difference(m, min, v_2, v_3);
	
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*min + t];
	
					//t = 0;
				}
	
	
				//row difference for last row
				i = last - 1;
	
				p_a = a + (i*n+s)*m;
	
				t = 0;
	
				for( s1 = 0; s1 < l; s1++)
				{
					for(q = 0; q < m; q++) 
					{
						v_1[t] = p_a[ s1*n + q ];
						p_a[ s1*n + q ] = 0.0;
						t++;
					}
				}
				//t = 0;
	
				for(s1 = 0; s1 < min; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
	
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
	
				//matrix_difference(min, 1, v_1, v_3);// V_2 = V_2 - V_3
	
				for(s1 = 0; s1 < min; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
	
				//t = 0;
	
				for( j = s + 1; j < last - 1; j ++ )
				{
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
					//p_a = a + (s*n+j)*m;
	
						//matrix_multiplication(min, m, m, v_1, v_2, v_3);
	
						p_a = a + (s_n + j)*m;
	
						for(s1 = 0; s1 < min; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
	
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
								v_2[s1*m + t] = sum;
							}
						}
						//t = 0;
	
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
						p_a = a + (i*n+j)*m;
	
						//matrix_difference(min, m, v_2, v_3);
	
						for(s1 = 0; s1 < min; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					
				}
	
				j = last - 1;
					
						//matrix_multiplication(min, m, min, v_1, v_2, v_3);
				p_a = a + (s_n + j)*m;
	
				for(s1 = 0; s1 < min; s1 ++)
				{
					for(t = 0; t < min; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
						v_2[s1*min + t] = sum;
					}
				}
	
				p_a = a + (i*n+j)*m;
	
				//matrix_difference(min, min, v_2, v_3);
	
				for(s1 = 0; s1 < min; s1 ++)
					for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*min + t];
	
				//t = 0;
				// END do row difference
			}
	
			// inverse for last row
			s = last - 1;
	
			s_n = s * n;
	
			p_a = a + (s_n +s )*m;
	
			last1 = l;
	
			for(q = 0; q < min; q ++ )
			{
				for( d = 0; d < min; d ++ ) 
				{
					if( d == q) v_1[q*min + d] = 1.0;
					else v_1[q*min + d] = 0.0;
				}
			}
	
			for(d = 0; d < min; d++ )  // step number s
			{
				//print_matrix(n, n, r, a);
				//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
				if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
				{
					printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
					printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
					return 1;
				}
	
				// futher a_ss != 0
				divisor = p_a[d*n + d];
				//printf("\ndivisor is %10.3e\n", divisor);
	
				for(t = 0; t < min; t++) // made a_ss equal to 1
				{
					p_a[d*n + t] = p_a[d*n + t] / divisor; 
	
					v_1[d*min + t] = v_1[d*min + t] / divisor; 
				}
	
				//print_matrix(n, n, r, v_1);
	
				for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
				{
					//if( s1 != d) 
					//{
						if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
						{
							divisor = p_a[s1*n + d];
	
							//printf("\ncycle divisor is %10.3e\n", divisor);
	
							for( t = 0; t < min; t++ )	
							{
								//printf("\n%10.10Le = %10.20Le - %10.20Le * %10.20Le\n", v_1[s1*min + t] - divisor * v_1[d*min + t], v_1[s1*min + t], divisor, v_1[d*min + t]);
	
								p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
	
								v_1[s1*min + t] = v_1[s1*min + t] - divisor * v_1[d*min + t];
							}
	
						}
					//}
				}
				for(s1 = d + 1; s1 < min; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
				{
					//if( s1 != d) 
					//{
						if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
						{
							divisor = p_a[s1*n + d];
	
							//printf("\ncycle divisor is %10.3e\n", divisor);
	
							for( t = 0; t < min; t++ )	
							{
								//printf("\n%10.10Le = %10.20Le - %10.20Le * %10.20Le\n", v_1[s1*min + t] - divisor * v_1[d*min + t], v_1[s1*min + t], divisor, v_1[d*min + t]);
	
								p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
	
								v_1[s1*min + t] = v_1[s1*min + t] - divisor * v_1[d*min + t];
							}
	
						}
					//}
				}
			}
	
			//t = 0;
	
			for(s1 = 0; s1 < min; s1 ++)
			{
				for(t = 0; t < 1; t ++ )
				{
					sum = 0.0;
	
					for(q = 0; q < min; q ++ ) sum += (v_1[s1*min + q]) * (b[s*m + q*1 + t]);
	
					v_2[s1*1 + t] = sum;
				}
		    }
			//t = 0;
	
			//put_back_i_vector_part(s, m, k, l, b, v_3); // V_3 and V_1 are free
	
			/*if( l != 0)
			{
				if( s == k) last1 = l;
				else last1 = m;
			}
			else last1 = m;*/
	
			for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
	
			//print_matrix(n, n, 7, a);
	
			//for(i = 0; i < 7; i++) printf(" %10.3e\n", b[i]);
	
			// BEGIN of X_i
		
			for( i = 0; i < min; i ++)	x[(last - 1)*m + i] = b[(last - 1)*m + i];
		
			for( i = last - 2; i >= 0; i --)
			{
		
				//return_i_vector_part(last - 1, m, k, l, b, v_1);
		
				p_a = a + (i*n+(last - 1))*m;
		
				// matrix_multiplication(m, m, 1, v_2, v_1, v_3);
				
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < min; q ++ ) sum += (p_a[s1*n + q]) * (b[(last - 1)*m + q*1 + t]);//(v_2[s1*m + q]) * (v_1[q*1 + t]);
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
				//matrix_difference(m, 1, v_1, v_3);
		
				for(s1 = 0; s1 < m; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
				//t = 0;
		
				for( j = last - 2; j > i ; j --)
				{ 
		
					//return_i_vector_part(j, m, k, l, x, v_1);// put X_i into V_2
		
					p_a = a + (i*n+j)*m;
		
					//matrix_multiplication(m, m, 1, v_2, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[j*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
					//t = 0;
		
					//matrix_difference(m, 1, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
				}
		
				for( j = 0; j < m; j ++)	x[i*m + j] = b[i*m + j];
		
			}
		}	
	}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////




















/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(m % 3 == 2)
	{
		m_new = m - 2;

		if(l == 0)
		{
			a_v = m;
			a_h = m;
	
			for(s = 0; s < last; s ++)
			{
				// BEGIN making A_ss E
	
				s_n = s*n;
		
				p_a = a + (s_n+s)*m;
		
				last1 = m;
		
				
					//return_i_j_matrix_block(s, s, n, m, k, l, a, v_1);// put A_ss into V_1
		
					//inverse_matrix(m, v_1, v_2); // V_2 = ( A_ss )^(-1)
		
				divisor = 0.0;
		
				for(q = 0; q < m; q ++ )
				{
					for( d = 0; d < m; d ++ ) 
					{
						if( d == q) v_1[q*m + d] = 1.0;
						else v_1[q*m + d] = 0.0;
					}
				}
		
				for(d = 0; d < m; d++ )  // step number s
				{	
					//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
					if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
					{
						printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
						printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
						return 1;
					}
		
					// futher a_ss != 0
					divisor = p_a[d*n + d];
		
					for(t = 0; t < m; t++) // made a_ss equal to 1
					{
						p_a[d*n + t] = p_a[d*n + t] / divisor; 
		
						v_1[d*m + t] = v_1[d*m + t] / divisor; 
					}
		
					for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
					for(s1 = d + 1; s1 < m; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
				}
		
				//t = 0;
		
				//matrix_multiplication(m, m, 1, v_2, v_1, v_3);// ( A_ss )^(-1) * B_s
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);;
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
					//put_back_i_vector_part(s, m, k, l, b, v_3); // V_2, V_3 and V_1 are free
		
					/*if( l != 0)
					{
						if( s == k) last1 = l;
						else last1 = m;
					}
					else last1 = m;*/
		
				for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
		
				// END making A_ss E
				// BEGIN ( A_ss )^(-1) * A_sj
		
				for( j = s + 1; j < last; j ++ )
				{
					p_a = a + (s_n+j)*m;
					
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_1);// put A_ij into V_1
		
					//matrix_multiplication(m, m, m, v_2, v_1, v_3);// V_3 = ( A_ss )^(-1) * A_sj
		
					/*for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < m; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]);
		
							v_2[s1*m + t] = sum;
						}
					}*/


					//print_matrix(9, 9, 9, v_2);


					for(r = 0; r < m_new; r += 3)
					{
						r_m = r * m;
				
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
				
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
				
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
				
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
								c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
							}
				
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
							
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
				
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
							v_2[r_m + 2*m + y + 2] = c22;
				
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
				
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
				
							c20 = 0.0;
							c21 = 0.0;
							c22 = 0.0;
						}
				
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
				
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
				
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
							c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
				
							c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
							c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
						}
				
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
							
						v_2[r_m + m + y]     = c10;
						v_2[r_m + m + y + 1] = c11;
				
						v_2[r_m + 2*m + y]     = c20;
						v_2[r_m + 2*m + y + 1] = c21;
				
						c00 = 0.0;
						c01 = 0.0;
				
						c10 = 0.0;
						c11 = 0.0;
				
						c20 = 0.0;
						c21 = 0.0;
					}
				
					//for(r = m_new; r < m; r ++)
					//{
					r_m = m_new * m;
					
					for(y = 0; y < m_new; y += 3)
					{
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
					
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
							c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
							c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
							c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
						}
					
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
						v_2[r_m + y + 2] = c02;
								
						v_2[r_m + m + y]     = c10;
						v_2[r_m + m + y + 1] = c11;
						v_2[r_m + m + y + 2] = c12;
					
						c00 = 0.0;
						c01 = 0.0;
						c02 = 0.0;
					
						c10 = 0.0;
						c11 = 0.0;
						c12 = 0.0;
					}
						//}
					for(z = 0; z < m; z ++)
					{
						z_m = z * n;
					
						c00 += v_1[r_m + z] * p_a[z_m + y];
						c01 += v_1[r_m + z] * p_a[z_m + y + 1];
					
						c10 += v_1[r_m + m + z] * p_a[z_m + y];
						c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
					}
					
					v_2[r_m + y]     = c00;
					v_2[r_m + y + 1] = c01;
								
					v_2[r_m + m + y]     = c10;
					v_2[r_m + m + y + 1] = c11;
					
					c00 = 0.0;
					c01 = 0.0;
					
					c10 = 0.0;
					c11 = 0.0;
					//print_matrix(9, 9, 9, v_2);










					//t = 0;
		
					//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3);// V_3 -> A_sj
		
					/*if( s == k) a_v = l;
					else a_v = m;
					if( j == k) a_h = l;
					else a_h = m;*/
		
					p_a = a + (s_n+j)*m;
	
					t = 0;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							p_a[ s1*n + q ] = v_2[t];
							t++;
						}
					}
					//t = 0;
					
				}
		
				// END ( A_ss )^(-1) * A_sj
				// BEGIN do row difference 
		
				for( i = s + 1; i < last; i ++ )
				{
		
					//return_i_j_matrix_block(i, s, n, m, k, l, a, v_1);// put A_is into V_1
		
					t = 0;
		
					p_a = a + (i*n+s)*m;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							v_1[t] = p_a[ s1*n + q ];
							p_a[ s1*n + q ] = 0.0;
							t++;
						}
					}
					//t = 0;
		
					//matrix_multiplication(m, m, 1, v_1, v_2, v_2); 
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
		
						//matrix_difference(m, 1, v_2, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
		
					for( j = s + 1; j < last; j ++ )
					{
						//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (s_n+j)*m;
		
						//matrix_multiplication(m, m, m, v_1, v_2, v_3);
		
						/*for(s1 = 0; s1 < m; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
		
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
		
								v_2[s1*m + t] = sum;
							}
						}*/
						//t = 0;

						for(r = 0; r < m_new; r += 3)
						{
							r_m = r * m;
					
							for(y = 0; y < m_new; y += 3)
							{
								for(z = 0; z < m; z ++)
								{
									z_m = z * n;
					
									c00 += v_1[r_m + z] * p_a[z_m + y];
									c01 += v_1[r_m + z] * p_a[z_m + y + 1];
									c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
									c10 += v_1[r_m + m + z] * p_a[z_m + y];
									c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
									c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
									c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
									c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
									c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
								}
					
								v_2[r_m + y]     = c00;
								v_2[r_m + y + 1] = c01;
								v_2[r_m + y + 2] = c02;
								
								v_2[r_m + m + y]     = c10;
								v_2[r_m + m + y + 1] = c11;
								v_2[r_m + m + y + 2] = c12;
					
								v_2[r_m + 2*m + y]     = c20;
								v_2[r_m + 2*m + y + 1] = c21;
								v_2[r_m + 2*m + y + 2] = c22;
					
								c00 = 0.0;
								c01 = 0.0;
								c02 = 0.0;
					
								c10 = 0.0;
								c11 = 0.0;
								c12 = 0.0;
					
								c20 = 0.0;
								c21 = 0.0;
								c22 = 0.0;
							}
					
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
					
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
					
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
					
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
							}
					
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
								
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
					
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
					
							c00 = 0.0;
							c01 = 0.0;
					
							c10 = 0.0;
							c11 = 0.0;
					
							c20 = 0.0;
							c21 = 0.0;
						}
					
						//for(r = m_new; r < m; r ++)
						//{
						r_m = m_new * m;
						
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
						
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
						
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
						
							}
						
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
									
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
						
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
						
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
						}
							//}
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
						
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
						
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
							c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
						}
						
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
									
						v_2[r_m + m + y]     = c10;
						v_2[r_m + m + y + 1] = c11;
						
						c00 = 0.0;
						c01 = 0.0;
						
						c10 = 0.0;
						c11 = 0.0;




















		
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (i*n+j)*m;
		
						//matrix_difference(m, m, v_2, v_3);
		
						for(s1 = 0; s1 < m; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					}
				}
				// END do row difference
			}
			
			// BEGIN of X_i
		
			for( i = 0; i < min; i ++)	x[(last - 1)*m + i] = b[(last - 1)*m + i];
		
			for( i = last - 2; i >= 0; i --)
			{
		
				//return_i_vector_part(last - 1, m, k, l, b, v_1);
		
				p_a = a + (i*n+(last - 1))*m;
		
				// matrix_multiplication(m, m, 1, v_2, v_1, v_3);
				
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[(last - 1)*m + q*1 + t]);//(v_2[s1*m + q]) * (v_1[q*1 + t]);
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
				//matrix_difference(m, 1, v_1, v_3);
		
				for(s1 = 0; s1 < m; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
				//t = 0;
		
				for( j = last - 2; j > i ; j --)
				{ 
		
					//return_i_vector_part(j, m, k, l, x, v_1);// put X_i into V_2
		
					p_a = a + (i*n+j)*m;
		
					//matrix_multiplication(m, m, 1, v_2, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[j*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
					//t = 0;
		
					//matrix_difference(m, 1, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
				}
		
				for( j = 0; j < m; j ++)	x[i*m + j] = b[i*m + j];
		
			}
		}
	
		else // l != 0
		{
			a_v = m;
			a_h = m;
	
			for(s = 0; s < last - 1; s ++)
			{
				// BEGIN making A_ss E
	
				s_n = s * n;
		
				p_a = a + (s_n+s)*m;
		
				last1 = m;
		
				
					//return_i_j_matrix_block(s, s, n, m, k, l, a, v_1);// put A_ss into V_1
		
					//inverse_matrix(m, v_1, v_2); // V_2 = ( A_ss )^(-1)
		
				for(q = 0; q < m; q ++ )
				{
					for( d = 0; d < m; d ++ ) 
					{
						if( d == q) v_1[q*m + d] = 1.0;
						else v_1[q*m + d] = 0.0;
					}
				}
		
				for(d = 0; d < m; d++ )  // step number s
				{	
					//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
					if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
					{
						printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
						printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
						return 1;
					}
		
					// futher a_ss != 0
					divisor = p_a[d*n + d];
		
					for(t = 0; t < m; t++) // made a_ss equal to 1
					{
						p_a[d*n + t] = p_a[d*n + t] / divisor; 
		
						v_1[d*m + t] = v_1[d*m + t] / divisor; 
					}
		
					for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
					for(s1 = d + 1; s1 < m; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
					{
						//if( s1 != d) 
						//{
							if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
							{
								divisor = p_a[s1*n + d];
		
								for( t = 0; t < m; t++ )	
								{
									p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
		
									v_1[s1*m + t] = v_1[s1*m + t] - divisor * v_1[d*m + t];
								}
		
							}
						//}
					}
				}
		
				//t = 0;
		
				//matrix_multiplication(m, m, 1, v_2, v_1, v_3);// ( A_ss )^(-1) * B_s
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);;
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
					//put_back_i_vector_part(s, m, k, l, b, v_3); // V_2, V_3 and V_1 are free
		
					/*if( l != 0)
					{
						if( s == k) last1 = l;
						else last1 = m;
					}
					else last1 = m;*/
		
				for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
		
				// END making A_ss E
				// BEGIN ( A_ss )^(-1) * A_sj
		
				for( j = s + 1; j < last - 1; j ++ )
				{
					p_a = a + (s_n+j)*m;
					
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_1);// put A_ij into V_1
		
					//matrix_multiplication(m, m, m, v_2, v_1, v_3);// V_3 = ( A_ss )^(-1) * A_sj
		
					/*for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < m; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]);
		
							v_2[s1*m + t] = sum;
						}
					}*/






					for(r = 0; r < m_new; r += 3)
					{
						r_m = r * m;
				
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
				
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
				
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
				
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
								c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
							}
				
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
							
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
				
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
							v_2[r_m + 2*m + y + 2] = c22;
				
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
				
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
				
							c20 = 0.0;
							c21 = 0.0;
							c22 = 0.0;
						}
				
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
				
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
				
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
							c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
				
							c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
							c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
						}
				
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
							
						v_2[r_m + m + y]     = c10;
						v_2[r_m + m + y + 1] = c11;
				
						v_2[r_m + 2*m + y]     = c20;
						v_2[r_m + 2*m + y + 1] = c21;
				
						c00 = 0.0;
						c01 = 0.0;
				
						c10 = 0.0;
						c11 = 0.0;
				
						c20 = 0.0;
						c21 = 0.0;
					}
				
					//for(r = m_new; r < m; r ++)
					//{
					r_m = m_new * m;
					
					for(y = 0; y < m_new; y += 3)
					{
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
					
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
							c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
							c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
							c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
						}
					
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
						v_2[r_m + y + 2] = c02;
								
						v_2[r_m + m + y]     = c10;
						v_2[r_m + m + y + 1] = c11;
						v_2[r_m + m + y + 2] = c12;
					
						c00 = 0.0;
						c01 = 0.0;
						c02 = 0.0;
					
						c10 = 0.0;
						c11 = 0.0;
						c12 = 0.0;
					}
						//}
					for(z = 0; z < m; z ++)
					{
						z_m = z * n;
					
						c00 += v_1[r_m + z] * p_a[z_m + y];
						c01 += v_1[r_m + z] * p_a[z_m + y + 1];
					
						c10 += v_1[r_m + m + z] * p_a[z_m + y];
						c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
					}
					
					v_2[r_m + y]     = c00;
					v_2[r_m + y + 1] = c01;
								
					v_2[r_m + m + y]     = c10;
					v_2[r_m + m + y + 1] = c11;
					
					c00 = 0.0;
					c01 = 0.0;
					
					c10 = 0.0;
					c11 = 0.0;
					//t = 0;
		
					//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3);// V_3 -> A_sj
		
					/*if( s == k) a_v = l;
					else a_v = m;
					if( j == k) a_h = l;
					else a_h = m;*/
		
					//p_a = a + (s*n+j)*m;
					t = 0;
		
					for( s1 = 0; s1 < m; s1++)
					{
						for(q = 0; q < m; q++) 
						{
							p_a[ s1*n + q ] = v_2[t];
							t++;
						}
					}
					//t = 0;
					
				}
	
				// j == last - 1
				j = last - 1;
	
				p_a = a + (s_n + j)*m;
	
	
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < min; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[ q*n + t ]); //(v_1[q*min + t]);
	
						v_2[s1*min + t] = sum;
					}
				}
				//t = 0;
	
				//put_back_i_j_matrix_block(s, j, n, m, k, l, a, v_3); // V_3 -> A_sj
	
				/*if( s == k) a_v = l;
				else a_v = m;
				if( j == k) a_h = l;
				else a_h = m;*/
	
				p_a = a + (s_n + j)*m;
	
				t = 0;
	
				for( s1 = 0; s1 < m; s1++)
				{
					for(q = 0; q < l; q++) 
					{
						p_a[ s1*n + q ] = v_2[t];
						t++;
					}
				}
				//t = 0;
	
	
	
	
	
		
				// END ( A_ss )^(-1) * A_sj
				// BEGIN do row difference 
	
				a_h = m;
		
				for( i = s + 1; i < last - 1; i ++ )
				{
	
					//print_matrix(n, n, 7, a);
		
					//return_i_j_matrix_block(i, s, n, m, k, l, a, v_1);// put A_is into V_1
		
					t = 0;
		
					p_a = a + (i*n+s)*m;
		
					for( s1 = 0; s1 < a_v; s1++)
					{
						for(q = 0; q < a_h; q++) 
						{
							v_1[t] = p_a[ s1*n + q ];
							p_a[ s1*n + q ] = 0.0;
							t++;
						}
					}
					//t = 0;
		
					//matrix_multiplication(m, m, 1, v_1, v_2, v_2); 
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
		
						//matrix_difference(m, 1, v_2, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
		
					for( j = s + 1; j < last - 1; j ++ )
					{
						//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (s_n + j)*m;
		
						//matrix_multiplication(m, m, m, v_1, v_2, v_3);
		
						/*for(s1 = 0; s1 < m; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
		
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
		
								v_2[s1*m + t] = sum;
							}
						}*/


						for(r = 0; r < m_new; r += 3)
						{
							r_m = r * m;
					
							for(y = 0; y < m_new; y += 3)
							{
								for(z = 0; z < m; z ++)
								{
									z_m = z * n;
					
									c00 += v_1[r_m + z] * p_a[z_m + y];
									c01 += v_1[r_m + z] * p_a[z_m + y + 1];
									c02 += v_1[r_m + z] * p_a[z_m + y + 2];
					
									c10 += v_1[r_m + m + z] * p_a[z_m + y];
									c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
									c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
					
									c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
									c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
									c22 += v_1[r_m + 2*m + z] * p_a[z_m + y + 2];
								}
					
								v_2[r_m + y]     = c00;
								v_2[r_m + y + 1] = c01;
								v_2[r_m + y + 2] = c02;
								
								v_2[r_m + m + y]     = c10;
								v_2[r_m + m + y + 1] = c11;
								v_2[r_m + m + y + 2] = c12;
					
								v_2[r_m + 2*m + y]     = c20;
								v_2[r_m + 2*m + y + 1] = c21;
								v_2[r_m + 2*m + y + 2] = c22;
					
								c00 = 0.0;
								c01 = 0.0;
								c02 = 0.0;
					
								c10 = 0.0;
								c11 = 0.0;
								c12 = 0.0;
					
								c20 = 0.0;
								c21 = 0.0;
								c22 = 0.0;
							}
					
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
					
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
					
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
					
								c20 += v_1[r_m + 2*m + z] * p_a[z_m + y];
								c21 += v_1[r_m + 2*m + z] * p_a[z_m + y + 1];
							}
					
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
								
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
					
							v_2[r_m + 2*m + y]     = c20;
							v_2[r_m + 2*m + y + 1] = c21;
					
							c00 = 0.0;
							c01 = 0.0;
					
							c10 = 0.0;
							c11 = 0.0;
					
							c20 = 0.0;
							c21 = 0.0;
						}
					
						//for(r = m_new; r < m; r ++)
						//{
						r_m = m_new * m;
						
						for(y = 0; y < m_new; y += 3)
						{
							for(z = 0; z < m; z ++)
							{
								z_m = z * n;
						
								c00 += v_1[r_m + z] * p_a[z_m + y];
								c01 += v_1[r_m + z] * p_a[z_m + y + 1];
								c02 += v_1[r_m + z] * p_a[z_m + y + 2];
						
								c10 += v_1[r_m + m + z] * p_a[z_m + y];
								c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
								c12 += v_1[r_m + m + z] * p_a[z_m + y + 2];
						
							}
						
							v_2[r_m + y]     = c00;
							v_2[r_m + y + 1] = c01;
							v_2[r_m + y + 2] = c02;
									
							v_2[r_m + m + y]     = c10;
							v_2[r_m + m + y + 1] = c11;
							v_2[r_m + m + y + 2] = c12;
						
							c00 = 0.0;
							c01 = 0.0;
							c02 = 0.0;
						
							c10 = 0.0;
							c11 = 0.0;
							c12 = 0.0;
						}
							//}
						for(z = 0; z < m; z ++)
						{
							z_m = z * n;
						
							c00 += v_1[r_m + z] * p_a[z_m + y];
							c01 += v_1[r_m + z] * p_a[z_m + y + 1];
						
							c10 += v_1[r_m + m + z] * p_a[z_m + y];
							c11 += v_1[r_m + m + z] * p_a[z_m + y + 1];
						}
						
						v_2[r_m + y]     = c00;
						v_2[r_m + y + 1] = c01;
									
						v_2[r_m + m + y]     = c10;
						v_2[r_m + m + y + 1] = c11;
						
						c00 = 0.0;
						c01 = 0.0;
						
						c10 = 0.0;
						c11 = 0.0;








						//t = 0;
		
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
		
						p_a = a + (i*n+j)*m;
		
						//matrix_difference(m, m, v_2, v_3);
		
						for(s1 = 0; s1 < m; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					}
	
					//last column for not last row
	
					j = last - 1;
	
					p_a = a + (s_n + j)*m;
	
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < min; t ++ )
						{
							sum = 0.0;
	
							for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
							v_2[s1*min + t] = sum;
						}
					}
					//t = 0;
	
					//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
					p_a = a + (i*n+j)*m;
	
					//matrix_difference(m, min, v_2, v_3);
	
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*min + t];
	
					//t = 0;
				}
	
	
				//row difference for last row
				i = last - 1;
	
				p_a = a + (i*n+s)*m;
	
				t = 0;
	
				for( s1 = 0; s1 < l; s1++)
				{
					for(q = 0; q < m; q++) 
					{
						v_1[t] = p_a[ s1*n + q ];
						p_a[ s1*n + q ] = 0.0;
						t++;
					}
				}
				//t = 0;
	
				for(s1 = 0; s1 < min; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (b[s*m + q*1 + t]);
	
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
	
				//matrix_difference(min, 1, v_1, v_3);// V_2 = V_2 - V_3
	
				for(s1 = 0; s1 < min; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
	
				//t = 0;
	
				for( j = s + 1; j < last - 1; j ++ )
				{
					//return_i_j_matrix_block(s, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
					//p_a = a + (s*n+j)*m;
	
						//matrix_multiplication(min, m, m, v_1, v_2, v_3);
	
						p_a = a + (s_n + j)*m;
	
						for(s1 = 0; s1 < min; s1 ++)
						{
							for(t = 0; t < m; t ++ )
							{
								sum = 0.0;
	
								for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
								v_2[s1*m + t] = sum;
							}
						}
						//t = 0;
	
						//return_i_j_matrix_block(i, j, n, m, k, l, a, v_2);// put A_sj into V_2
	
						p_a = a + (i*n+j)*m;
	
						//matrix_difference(min, m, v_2, v_3);
	
						for(s1 = 0; s1 < min; s1 ++)
							for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
	
						//t = 0;
					
				}
	
				j = last - 1;
					
						//matrix_multiplication(min, m, min, v_1, v_2, v_3);
				p_a = a + (s_n + j)*m;
	
				for(s1 = 0; s1 < min; s1 ++)
				{
					for(t = 0; t < min; t ++ )
					{
						sum = 0.0;
	
						for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (p_a[q*n + t]);
	
						v_2[s1*min + t] = sum;
					}
				}
	
				p_a = a + (i*n+j)*m;
	
				//matrix_difference(min, min, v_2, v_3);
	
				for(s1 = 0; s1 < min; s1 ++)
					for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*min + t];
	
				//t = 0;
				// END do row difference
			}
	
			// inverse for last row
			s = last - 1;
	
			s_n = s * n;
	
			p_a = a + (s_n +s )*m;
	
			last1 = l;
	
			for(q = 0; q < min; q ++ )
			{
				for( d = 0; d < min; d ++ ) 
				{
					if( d == q) v_1[q*min + d] = 1.0;
					else v_1[q*min + d] = 0.0;
				}
			}
	
			for(d = 0; d < min; d++ )  // step number s
			{
				//print_matrix(n, n, r, a);
				//printf("\n%d*min + %d , %10.3Le <= %10.3Le\n", d, d, fabs(p_a[d*n + d]), eps * norm );
				if(fabs(p_a[d*n + d]) <= eps * norm)   // if a_ss == 0
				{
					printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
					printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(p_a[d*n + d]), eps * norm );
					return 1;
				}
	
				// futher a_ss != 0
				divisor = p_a[d*n + d];
				//printf("\ndivisor is %10.3e\n", divisor);
	
				for(t = 0; t < min; t++) // made a_ss equal to 1
				{
					p_a[d*n + t] = p_a[d*n + t] / divisor; 
	
					v_1[d*min + t] = v_1[d*min + t] / divisor; 
				}
	
				//print_matrix(n, n, r, v_1);
	
				for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
				{
					//if( s1 != d) 
					//{
						if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
						{
							divisor = p_a[s1*n + d];
	
							//printf("\ncycle divisor is %10.3e\n", divisor);
	
							for( t = 0; t < min; t++ )	
							{
								//printf("\n%10.10Le = %10.20Le - %10.20Le * %10.20Le\n", v_1[s1*min + t] - divisor * v_1[d*min + t], v_1[s1*min + t], divisor, v_1[d*min + t]);
	
								p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
	
								v_1[s1*min + t] = v_1[s1*min + t] - divisor * v_1[d*min + t];
							}
	
						}
					//}
				}
				for(s1 = d + 1; s1 < min; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
				{
					//if( s1 != d) 
					//{
						if( fabs(p_a[s1*n + d]) > eps * norm )  // if a_is != 0
						{
							divisor = p_a[s1*n + d];
	
							//printf("\ncycle divisor is %10.3e\n", divisor);
	
							for( t = 0; t < min; t++ )	
							{
								//printf("\n%10.10Le = %10.20Le - %10.20Le * %10.20Le\n", v_1[s1*min + t] - divisor * v_1[d*min + t], v_1[s1*min + t], divisor, v_1[d*min + t]);
	
								p_a[s1*n + t] = p_a[s1*n + t] - divisor * p_a[d*n + t];
	
								v_1[s1*min + t] = v_1[s1*min + t] - divisor * v_1[d*min + t];
							}
	
						}
					//}
				}
			}
	
			//t = 0;
	
			for(s1 = 0; s1 < min; s1 ++)
			{
				for(t = 0; t < 1; t ++ )
				{
					sum = 0.0;
	
					for(q = 0; q < min; q ++ ) sum += (v_1[s1*min + q]) * (b[s*m + q*1 + t]);
	
					v_2[s1*1 + t] = sum;
				}
		    }
			//t = 0;
	
			//put_back_i_vector_part(s, m, k, l, b, v_3); // V_3 and V_1 are free
	
			/*if( l != 0)
			{
				if( s == k) last1 = l;
				else last1 = m;
			}
			else last1 = m;*/
	
			for(s1 = 0; s1 < last1; s1++) b[s*m + s1] = v_2[s1];
	
			//print_matrix(n, n, 7, a);
	
			//for(i = 0; i < 7; i++) printf(" %10.3e\n", b[i]);
	
			// BEGIN of X_i
		
			for( i = 0; i < min; i ++)	x[(last - 1)*m + i] = b[(last - 1)*m + i];
		
			for( i = last - 2; i >= 0; i --)
			{
		
				//return_i_vector_part(last - 1, m, k, l, b, v_1);
		
				p_a = a + (i*n+(last - 1))*m;
		
				// matrix_multiplication(m, m, 1, v_2, v_1, v_3);
				
				for(s1 = 0; s1 < m; s1 ++)
				{
					for(t = 0; t < 1; t ++ )
					{
						sum = 0.0;
		
						for(q = 0; q < min; q ++ ) sum += (p_a[s1*n + q]) * (b[(last - 1)*m + q*1 + t]);//(v_2[s1*m + q]) * (v_1[q*1 + t]);
		
						v_2[s1*1 + t] = sum;
					}
				}
				//t = 0;
		
				//matrix_difference(m, 1, v_1, v_3);
		
				for(s1 = 0; s1 < m; s1 ++)
					for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
				//t = 0;
		
				for( j = last - 2; j > i ; j --)
				{ 
		
					//return_i_vector_part(j, m, k, l, x, v_1);// put X_i into V_2
		
					p_a = a + (i*n+j)*m;
		
					//matrix_multiplication(m, m, 1, v_2, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
					{
						for(t = 0; t < 1; t ++ )
						{
							sum = 0.0;
		
							for(q = 0; q < m; q ++ ) sum += (p_a[s1*n + q]) * (b[j*m + q*1 + t]);
		
							v_2[s1*1 + t] = sum;
						}
					}
					//t = 0;
		
					//matrix_difference(m, 1, v_1, v_3);
		
					for(s1 = 0; s1 < m; s1 ++)
						for(t = 0; t < 1; t ++ ) b[i*m + s1*1 + t] = b[i*m + s1*1 + t] - v_2[s1*1 + t];
		
					//t = 0;
				}
		
				for( j = 0; j < m; j ++)	x[i*m + j] = b[i*m + j];
		
			}
		}	
	}


	// END of X_i

	//gettimeofday(&tv2, &tz);
	//std::cout<<" Time of searching of the solution = "<<tv2.tv_sec+tv2.tv_usec/1000000.0 - tv1.tv_sec-tv1.tv_usec/1000000.0<<" sec"<<endl<<endl;

	return 0;
}


