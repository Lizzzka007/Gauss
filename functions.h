#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void vector(int n, double * a, double * b);

//void return_i_j_matrix_block(int i, int j, int n, int m, int k, int l, double *a, double * v_1);

//void put_back_i_j_matrix_block(int i, int j, int n, int m, int k, int l, double *a, double * v_1);// put V_1 into A_ij 

int print_matrix(int l, int n, int r, double * a);

void enter_matrix(int n, int m, int k, int s, int l, double * a);

int enter_matrix_from_file(double * a, int n, const char * filename);

//void inverse_matrix(int m, double * a, double * unit);

//void matrix_multiplication(int l, int m, int n, double * a, double * b, double * v);

//void null_matrix(int l, int m, double * a);

//void matrix_difference(int l, int m, double * a, double * b);

int solution( int n, int m, double *a, double *b, double *x, double *v_1, double *v_2);

//void return_i_vector_part(int i, int m, int k, int l, double *b, double * v_1);// put B_i into V_1

//void put_back_i_vector_part(int i, int m, int k, int l, double *b, double * v_1);// put V_1 into B_i

void print_discrepancy(int s5, int n,int m, double * a, double * b, double * x, double *v_1, double *v_2, float tv_seconds, const char * name);








#endif
